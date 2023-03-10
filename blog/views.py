from django.views import generic
from .models import Posts, Category, Comments
from hitcount.views import HitCountDetailView
from django.views.generic.dates import MonthArchiveView
from django.shortcuts import get_object_or_404, reverse, redirect
from django.conf import settings
from django.db.models import Q
from .forms import AddCommentForm, ContactForm, ProfileModelForm
from django.http import HttpResponseRedirect
from django.contrib import messages
from accounts.models import Profile, ContactModel


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
        messages.success(request, "Yorumunuz ba??ar??l?? bir ??ekilde eklendi")
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


class PostDeleteView(generic.DeleteView):
    template_name = "blog/pages/delete.html"
    model = Posts

    def get(self, request, *args, **kwargs):
        if not self.request.user.is_staff:
            messages.error(request, "Yetkili giri??i yap??n??z!")
            return redirect('%s?next=/blog/' % (settings.LOGIN_URL))
        return super().get(request, *args, **kwargs)

    def get_success_url(self):
        messages.success(self.request, "G??nderi ba??ar??l?? bir ??ekilde silindi")
        return reverse('blog:anasayfa')


class CategoryDeleteView(generic.DeleteView):
    template_name = "blog/pages/delete.html"
    model = Category

    def get(self, request, *args, **kwargs):
        if not request.user.is_staff:
            messages.error(request, "Yetkili giri??i yap??n??z!")
            return redirect('%s?next=/blog/' % (settings.LOGIN_URL))
        return super().get(request, *args, **kwargs)

    def get_success_url(self):
        messages.success(self.request, "G??nderi ba??ar??l?? bir ??ekilde silindi")
        return reverse('blog:anasayfa')


class CommentDetailView(generic.DetailView):
    template_name = "blog/pages/comments.html"
    model = Comments

    def get(self, request, *args, **kwargs):
        if not request.user.is_staff:
            messages.error(request, "Yetkili giri??i yap??n??z!")
            return redirect('%s?next=/blog/' % (settings.LOGIN_URL))
        return super().get(request, *args, **kwargs)


def comment_read(request, pk):
    if request.user.is_anonymous:
        messages.error(request, "G??nderiyi be??enmek i??in giri?? yap??n??z!")
        return redirect('%s?next=/blog/' % (settings.LOGIN_URL))
    if not request.user.is_staff:
        messages.error(request, "Yetkili giri??i yap??n??z!")
        return redirect('%s?next=/blog/' % (settings.LOGIN_URL))
    comment = get_object_or_404(Comments, pk=pk)
    comment.status = "okundu"
    comment.save()
    return HttpResponseRedirect(reverse("dashboard:anasayfa", args=[request.user.username]))


def comment_delete(request, pk):
    if request.user.is_anonymous:
        messages.error(request, "G??nderiyi be??enmek i??in giri?? yap??n??z!")
        return redirect('%s?next=/blog/' % (settings.LOGIN_URL))
    comment = get_object_or_404(Comments, pk=pk)
    if not request.user.username == comment.commentator.username:
        messages.error(request, "Yetkisiz kullan??c??")
        return redirect('%s?next=/blog/' % (settings.LOGIN_URL))
    comment.delete()
    messages.success(request, "Yorumunuz ba??ar??l?? bir ??ekilde silindi")
    return HttpResponseRedirect(request.META.get('HTTP_REFERER'))


def comment_like(request, id):
    if request.user.is_anonymous:
        messages.error(request, "G??nderiyi be??enmek i??in giri?? yap??n??z!")
        return redirect('%s?next=/blog/' % (settings.LOGIN_URL))
    comment = get_object_or_404(Comments, id=id)
    if comment.commentator.username == request.user.username:
        messages.info(request, "Kendini be??enmi?? birisin :) ")
    else:
        messages.success(request, f"{comment.commentator.username} isimli kullan??c??n??n yorumunu be??endiniz")

    comment.likes.add(request.user)

    return HttpResponseRedirect(request.META.get('HTTP_REFERER'))


def comment_dislike(request, id):
    if request.user.is_anonymous:
        messages.error(request, "G??nderiyi be??enmek i??in giri?? yap??n??z!")
        return redirect('%s?next=/blog/' % (settings.LOGIN_URL))
    comment = get_object_or_404(Comments, id=id)
    comment.dislike.add(request.user)
    messages.error(request, f"{comment.commentator.username} isimli kullan??c??n??n yorumunu be??enmediniz")
    return HttpResponseRedirect(request.META.get('HTTP_REFERER'))


def report_comment(request, id):
    if request.user.is_anonymous:
        messages.error(request, "G??nderiyi be??enmek i??in giri?? yap??n??z!")
        return redirect('%s?next=/blog/' % (settings.LOGIN_URL))
    comment = get_object_or_404(Comments, id=id)
    comment.report.add(request.user)
    if comment.report.count() > 5:
        comment.delete()
    messages.success(request, "Geri bildiriminiz iletilmi??tir")
    return HttpResponseRedirect(request.META.get('HTTP_REFERER'))


def like_post(request, pk, slug):
    if request.user.is_anonymous:
        messages.error(request, "G??nderiyi be??enmek i??in giri?? yap??n??z!")
        return redirect('%s?next=/blog/' % (settings.LOGIN_URL))
    post = get_object_or_404(Posts, pk=pk, slug=slug)
    post.likes.add(request.user)
    messages.success(request, f"{post.title} ba??l??kl?? g??nderiyi be??endiniz")
    return HttpResponseRedirect(request.META.get('HTTP_REFERER'))


def dislike_post(request, pk, title):
    if request.user.is_anonymous:
        messages.error(request, "G??nderiyi be??enmek i??in giri?? yap??n??z!")
        return redirect('%s?next=/blog/' % (settings.LOGIN_URL))
    post = get_object_or_404(Posts, pk=pk, title=title)
    post.dislike.add(request.user)
    messages.error(request, f"{post.title} ba??l??kl?? g??nderiyi be??enmediniz")
    return HttpResponseRedirect(request.META.get('HTTP_REFERER'))


class CommentUpdate(generic.UpdateView):
    template_name = "blog/pages/update.html"
    model = Comments
    form_class = AddCommentForm

    def get_success_url(self):
        messages.success(self.request, "Yorumunuz ba??ar??l?? bir ??ekilde g??ncellendi..")
        return HttpResponseRedirect(
            self.request.build_absolute_uri()
        ).url


class ProfileView(generic.DetailView, generic.CreateView):
    http_method_names = ['post', 'get']
    template_name = 'blog/pages/profile.html'
    model = Profile
    form_class = ContactForm

    def post(self, request, *args, **kwargs):
        message = ContactModel(
            author=self.get_object(),
            sender=self.request.user,
            content=self.request.POST.get('content'),
            contact_email=self.request.POST.get('contact_email')
        )
        message.save()
        messages.success(request, "Mesaj??n??z yazara iletildi")
        return HttpResponseRedirect(self.request.build_absolute_uri())

    def get_context_data(self, *, object_list=None, **kwargs):
        context = super().get_context_data(**kwargs)
        for profile in Profile.objects.filter(user=self.kwargs['pk']):
            context['user_social'] = profile.user_social.values()
        return context


class ProfileUpdateViewNonStaff(generic.UpdateView):
    template_name = 'blog/pages/update.html'
    model = Profile
    form_class = ProfileModelForm

    def form_valid(self, form):
        instance = form.save(commit=False)
        instance.user = self.request.user
        instance.save()
        messages.success(self.request, "Profil g??ncellendi!")
        return HttpResponseRedirect(reverse('blog:profile', kwargs={
            'user': self.request.user, 'pk': self.request.user.pk
        }))


