from django.views import generic
from .models import Posts, Category, Comments
from hitcount.views import HitCountDetailView
from django.views.generic.dates import MonthArchiveView
from django.shortcuts import get_object_or_404, reverse, redirect
from django.conf import settings
from django.db.models import Q
from .forms.comment import AddCommentForm
from .forms.post import CreatePostModelForm, CreateCategoryModelForm
from django.http import HttpResponseRedirect
from django.contrib import messages


class CategoriesView(generic.ListView):
    template_name = 'blog/pages/category.html'
    model = Category
    paginate_by = 10
    context_object_name = "category_list"

    def get(self, request, *args, **kwargs):
        return super().get(request, *args, **kwargs)

    def get_queryset(self):
        search = self.request.GET.get('search', False)
        if search:
            return Category.objects.filter(Q(title__icontains=search))
        else:
            return Category.objects.all()


class CategoryView(generic.ListView):
    template_name = 'blog/pages/category.html'
    model = Posts
    paginate_by = 3
    context_object_name = "category_post"

    def get_queryset(self):
        category = get_object_or_404(Category, slug=self.kwargs['slug'])
        return category.post.all().order_by('-id')

    def get_context_data(self, *, object_list=None, **kwargs):
        context = super().get_context_data(**kwargs)
        context['category_name'] = get_object_or_404(Category, slug=self.kwargs['slug'])
        return context


class BlogHomeView(generic.ListView):
    template_name = "blog/pages/home.html"
    model = Posts
    paginate_by = 10
    ordering = "-id"

    def get_context_data(self, *, object_list=None, **kwargs):
        context = super().get_context_data(**kwargs)
        if Posts.objects.exists():
            context['latest'] = Posts.objects.all().latest('created')
        context['popular_posts'] = Posts.objects.order_by('-hit_count__hits')[:2]
        context['popular_post_side'] = Posts.objects.order_by('-hit_count__hits')[2:10]
        context['archives'] = Posts.objects.dates('created', 'month', 'DESC')
        return context


class PostDetailView(HitCountDetailView, generic.DetailView):
    template_name = "blog/pages/detail.html"
    model = Posts
    count_hit = True
    form_class = AddCommentForm

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['popular_posts'] = Posts.objects.order_by('-hit_count__hits')[2:10]
        context['archives'] = Posts.objects.dates('created', 'month', 'DESC')
        posts = get_object_or_404(Posts, slug=self.kwargs['slug'])
        context['comments'] = posts.comment.all().order_by('-id')

        if self.request.user.is_authenticated:
            context['form'] = AddCommentForm(instance=self.request.user)
            context['post_list'] = Posts.objects.all().filter(author=self.request.user)
        return context

    def post(self, request, *args, **kwargs):
        comment = Comments(comment=request.POST.get('comment'),
                           commentator=self.request.user,
                           post=self.get_object())
        comment.save()
        messages.success(request, "Yorumunuz başarılı bir şekilde eklendi")
        return HttpResponseRedirect(self.request.build_absolute_uri())


class ArchiveView(MonthArchiveView):
    template_name = "blog/pages/category.html"
    model = Posts
    date_field = "created"
    allow_future = True
    paginate_by = 10
    month_format = "%m"

    def get_queryset(self):
        return Posts.objects.order_by('-id')

    def get_context_data(self, *, object_list=None, **kwargs):
        context = super().get_context_data(**kwargs)
        context['month'] = self.get_month()
        context['year'] = self.get_year()
        return context


class CreatePost(generic.CreateView):
    template_name = "dashboard/pages/create.html"
    model = Posts
    form_class = CreatePostModelForm

    def get(self, request, *args, **kwargs):
        if not self.request.user.is_staff:
            messages.error(request, "Yetkili girişi yapınız!")
            return redirect('%s?next=/blog/' % (settings.LOGIN_URL))
        return super().get(request, *args, **kwargs)

    def form_valid(self, form):
        instance = form.save(commit=False)
        instance.author = self.request.user
        instance.save()
        return super().form_valid(form)

    def get_success_url(self):
        messages.success(self.request, "Gönderi Başarılı bir şekilde yayınlandı..")
        return reverse('blog:anasayfa')


class CreateCategoryView(generic.CreateView):
    template_name = "dashboard/pages/create.html"
    model = Category
    form_class = CreateCategoryModelForm

    def get(self, request, *args, **kwargs):
        if not self.request.user.is_staff:
            messages.error(request, "Yetkili girişi yapınız!")
            return redirect('%s?next=/blog/' % (settings.LOGIN_URL))
        return super().get(request, *args, **kwargs)

    def get_success_url(self, *args, **kwargs):
        messages.success(self.request, "Kategori başarılı bir şekilde eklendi..")
        return HttpResponseRedirect(reverse("blog:post_create")).url


class PostUpdateView(generic.UpdateView):
    template_name = "dashboard/pages/update.html"
    model = Posts
    form_class = CreatePostModelForm

    def get(self, request, *args, **kwargs):
        if not self.request.user.is_staff:
            messages.error(request, "Yetkili girişi yapınız!")
            return redirect('%s?next=/blog/' % (settings.LOGIN_URL))
        return super().get(request, *args, **kwargs)

    def form_valid(self, form):
        instance = form.save(commit=False)
        instance.author = self.request.user
        instance.save()
        return super().form_valid(form)

    def get_success_url(self):
        messages.success(self.request, "Gönderi Başarılı bir şekilde güncellendi..")
        return reverse('blog:anasayfa')


class CategoryUpdateView(generic.UpdateView):
    template_name = "dashboard/pages/update.html"
    model = Category
    form_class = CreateCategoryModelForm

    def get(self, request, *args, **kwargs):
        if not self.request.user.is_staff:
            messages.error(request, "Yetkili girişi yapınız!")
            return redirect('%s?next=/blog/' % (settings.LOGIN_URL))
        return super().get(request, *args, **kwargs)

    def form_valid(self, form):
        instance = form.save(commit=False)
        instance.author = self.request.user
        instance.save()
        return super().form_valid(form)

    def get_success_url(self):
        messages.success(self.request, "Kategori başarılı bir şekilde güncellendi..")
        return reverse('blog:categories')


class PostDeleteView(generic.DeleteView):
    template_name = "dashboard/pages/delete.html"
    model = Posts

    def get(self, request, *args, **kwargs):
        if not self.request.user.is_staff:
            messages.error(request, "Yetkili girişi yapınız!")
            return redirect('%s?next=/blog/' % (settings.LOGIN_URL))
        return super().get(request, *args, **kwargs)

    def get_success_url(self):
        messages.success(self.request, "Gönderi başarılı bir şekilde silindi")
        return reverse('blog:anasayfa')


class CategoryDeleteView(generic.DeleteView):
    template_name = "dashboard/pages/delete.html"
    model = Category

    def get(self, request, *args, **kwargs):
        if not request.user.is_staff:
            messages.error(request, "Yetkili girişi yapınız!")
            return redirect('%s?next=/blog/' % (settings.LOGIN_URL))
        return super().get(request, *args, **kwargs)

    def get_success_url(self):
        messages.success(self.request, "Gönderi başarılı bir şekilde silindi")
        return reverse('blog:anasayfa')


class CommentDetailView(generic.DetailView):
    template_name = "dashboard/pages/comments.html"
    model = Comments

    def get(self, request, *args, **kwargs):
        if not request.user.is_staff:
            messages.error(request, "Yetkili girişi yapınız!")
            return redirect('%s?next=/blog/' % (settings.LOGIN_URL))
        return super().get(request, *args, **kwargs)


def comment_read(request, pk):
    if not request.user.is_staff:
        messages.error(request, "Yetkili girişi yapınız!")
        return redirect('%s?next=/blog/' % (settings.LOGIN_URL))
    comment = get_object_or_404(Comments, pk=pk)
    comment.status = "okundu"
    comment.save()
    return HttpResponseRedirect(reverse("dashboard:anasayfa", args=[request.user.username]))


def comment_delete(request, pk):
    comment = get_object_or_404(Comments, pk=pk)
    if not request.user.username == comment.commentator.username:
        messages.error(request, "Yetkisiz kullanıcı")
        return redirect('%s?next=/blog/' % (settings.LOGIN_URL))
    comment.delete()
    messages.success(request, "Yorumunuz başarılı bir şekilde silindi")
    return HttpResponseRedirect(reverse("blog:post_detail", args=(
        comment.post.category.title, comment.post.slug, comment.post.pk, comment.post.author, comment.post.author.pk,
        comment.post.created.date()
    )))


def comment_like(request, id):
    comment = get_object_or_404(Comments, id=id)
    if comment.commentator.username == request.user.username:
        messages.info(request, "Kendini beğenmiş birisin :) ")
    else:
        messages.success(request, f"{comment.commentator.username} isimli kullanıcının yorumunu beğendiniz")

    comment.likes.add(request.user)

    return HttpResponseRedirect(request.META.get('HTTP_REFERER'))


def comment_dislike(request, id):
    comment = get_object_or_404(Comments, id=id)
    comment.dislike.add(request.user)
    messages.error(request, f"{comment.commentator.username} isimli kullanıcının yorumunu beğenmediniz")
    return HttpResponseRedirect(reverse("blog:post_detail", args=(
        comment.post.category.title, comment.post.slug, comment.post.pk, comment.post.author, comment.post.author.pk,
        comment.post.created.date()
    )))


def report_comment(request, id):
    comment = get_object_or_404(Comments, id=id)
    comment.report.add(request.user)
    if comment.report.count() > 5:
        comment.delete()
    messages.success(request, "Geri bildiriminiz iletilmiştir")
    return HttpResponseRedirect(reverse("blog:post_detail", args=(
        comment.post.category.title, comment.post.slug, comment.post.pk, comment.post.author, comment.post.author.pk,
        comment.post.created.date()
    )))


def like_post(request, pk, title):
    if request.user.is_anonymous:
        messages.error(request, "Gönderiyi beğenmek için giriş yapınız!")
        return redirect('%s?next=/blog/' % (settings.LOGIN_URL))
    post = get_object_or_404(Posts, pk=pk, title=title)
    post.likes.add(request.user)
    messages.success(request, f"{post.title} başlıklı gönderiyi beğendiniz")
    return HttpResponseRedirect(reverse("blog:post_detail", args=(
        post.category.title, post.slug, post.pk, post.author, post.author.pk,
        post.created.date()
    )))


def dislike_post(request, pk, title):
    post = get_object_or_404(Posts, pk=pk, title=title)
    post.dislike.add(request.user)
    messages.error(request, f"{post.title} başlıklı gönderiyi beğenmediniz")
    return HttpResponseRedirect(reverse("blog:post_detail", args=(
        post.category.title, post.slug, post.pk, post.author, post.author.pk,
        post.created.date()
    )))
