from django.views import generic
from .models import Posts, Category, Comments, Subscribe, BlogContactModel
from hitcount.views import HitCountDetailView
from django.views.generic.dates import MonthArchiveView
from django.shortcuts import get_object_or_404, reverse, redirect, render
from django.conf import settings
from django.db.models import Q
from .forms import AddCommentForm, ContactForm, ProfileModelForm, ContactProfileForm, SubscribeModelForm
from django.http import HttpResponseRedirect
from django.contrib import messages
from accounts.models import Profile, ContactModel
from django.contrib.auth.models import User
from django.core.paginator import Paginator
from blog.forms.contact import BlogContactForm
from django.core.mail import send_mail
from dashboard.models import Notifications


class CategoriesView(generic.ListView, generic.FormView):
    template_name = 'blog/pages/category.html'
    model = Category
    paginate_by = 10
    context_object_name = "category_list"
    form_class = SubscribeModelForm

    def post(self, request, *args, **kwargs):
        email = request.POST.get('email')

        if email:

            if Subscribe.objects.exists():
                for abone in Subscribe.objects.all():
                    if email in abone.email:
                        messages.error(request, "Bültene zaten abonesiniz")
                        return HttpResponseRedirect(self.request.build_absolute_uri())

                    else:
                        Subscribe.objects.create(email=email)

            else:
                Subscribe.objects.create(email=email)

        messages.success(request, "Başarılı bir şekilde abone oldunuz")

        return HttpResponseRedirect(self.request.build_absolute_uri())

    def get_queryset(self):
        search = self.request.GET.get('search', False)
        if search:
            return Category.objects.filter(Q(title__icontains=search))
        else:
            return Category.objects.all()


class CategoryView(generic.ListView, generic.FormView):
    template_name = 'blog/pages/category.html'
    model = Posts
    paginate_by = 3
    context_object_name = "category_post"
    form_class = SubscribeModelForm

    def post(self, request, *args, **kwargs):

        email = request.POST.get('email')

        if Subscribe.objects.exists():
            for abone in Subscribe.objects.all():
                if email in abone.email:
                    messages.error(request, "Bültene zaten abonesiniz")
                    return HttpResponseRedirect(self.request.build_absolute_uri())

                else:
                    Subscribe.objects.create(email=email)

        else:
            Subscribe.objects.create(email=email)

        messages.success(request, "Başarılı bir şekilde abone oldunuz")

        return HttpResponseRedirect(self.request.build_absolute_uri())

    def get_queryset(self):
        category = get_object_or_404(Category, slug=self.kwargs['slug'])
        return category.post.all().order_by('-id')

    def get_context_data(self, *, object_list=None, **kwargs):
        context = super().get_context_data(**kwargs)
        context['category_name'] = get_object_or_404(Category, slug=self.kwargs['slug'])
        return context


class BlogHomeView(generic.ListView, generic.FormView):
    template_name = "blog/pages/home.html"
    model = Posts
    paginate_by = 10
    ordering = "-id"
    form_class = SubscribeModelForm

    def post(self, request, *args, **kwargs):
        email = request.POST.get('email')

        if Subscribe.objects.exists():
            for abone in Subscribe.objects.all():
                if email in abone.email:
                    messages.error(request, "Bültene zaten abonesiniz")
                    return HttpResponseRedirect(self.request.build_absolute_uri())

                else:
                    Subscribe.objects.create(email=email)

        else:
            Subscribe.objects.create(email=email)

        messages.success(request, "Başarılı bir şekilde abone oldunuz")

        return HttpResponseRedirect(self.request.build_absolute_uri())

    def get_context_data(self, *, object_list=None, **kwargs):
        context = super().get_context_data(**kwargs)
        if Posts.objects.exists():
            context['latest'] = Posts.objects.all().latest('created')
        context['popular_posts'] = Posts.objects.order_by('-hit_count__hits')[:2]
        context['popular_post_side'] = Posts.objects.order_by('-hit_count__hits')[2:10]
        context['archives'] = Posts.objects.dates('created', 'month', 'DESC')
        page = context['page_obj']
        paginator = page.paginator
        pagelist = paginator.get_elided_page_range(page.number, on_each_side=2, on_ends=2)
        context['pagelist'] = pagelist
        return context


class PostDetailView(HitCountDetailView, generic.DetailView, generic.FormView):
    template_name = "blog/pages/detail.html"
    model = Posts
    count_hit = True
    form_class = AddCommentForm

    def post(self, request, *args, **kwargs):
        email = request.POST.get('email')
        comment = request.POST.get('comment')

        if email:
            for abone in Subscribe.objects.all():
                if request.POST.get('email') in abone.email:
                    messages.success(request, "Bültene zaten abonesiniz")
                    return HttpResponseRedirect(self.request.build_absolute_uri())
                else:
                    Subscribe(email=email).save()
            messages.success(request, "Başarılı bir şekilde abone oldunuz")

        if comment:
            comment = Comments(comment=request.POST.get('comment'),
                               commentator=self.request.user,
                               post=self.get_object())
            send_mail(
                "Gönderinize yorum yapıldı",
                f"{self.request.user} isimli kullanıcı {comment.post.title} başlıklı gönderinize yorum yaptı.\n\n\n"
                f"Yorum: \n\n {comment.comment}."
                f"\n\n {redirect(self.request.build_absolute_uri()).url}",
                "bioblogdestek@gmail.com",
                [comment.post.author.email],
                fail_silently=False,
            )
            comment.save()

            messages.success(request, "Yorumunuz başarılı bir şekilde eklendi")

            if not Notifications.objects.filter(user=comment.post.author, type="yorum",
                                         title=f"{comment.commentator} {comment.post.title} başlıklı gönderinize yorum yaptı",
                                         text=comment.comment,
                                         url=redirect(self.request.build_absolute_uri()).url
                                         ):

                Notifications.objects.create(user=comment.post.author, type="yorum",
                                             title=f"{comment.commentator} {comment.post.title} başlıklı gönderinize yorum yaptı",
                                             text=comment.comment,
                                             url=redirect(self.request.build_absolute_uri()).url
                                             )

        return HttpResponseRedirect(self.request.build_absolute_uri())

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


class ArchiveView(MonthArchiveView, generic.FormView):
    template_name = "blog/pages/category.html"
    model = Posts
    date_field = "created"
    allow_future = True
    paginate_by = 10
    month_format = "%m"
    form_class = SubscribeModelForm

    def post(self, request, *args, **kwargs):
        email = request.POST.get('email')

        if Subscribe.objects.exists():
            for abone in Subscribe.objects.all():
                if email in abone.email:
                    messages.error(request, "Bültene zaten abonesiniz")
                    return HttpResponseRedirect(self.request.build_absolute_uri())

                else:
                    Subscribe.objects.create(email=email)

        else:
            Subscribe.objects.create(email=email)

        messages.success(request, "Başarılı bir şekilde abone oldunuz")

        return HttpResponseRedirect(self.request.build_absolute_uri())

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
            messages.error(request, "Yetkili girişi yapınız!")
            return redirect('%s?next=/blog/' % (settings.LOGIN_URL))
        return super().get(request, *args, **kwargs)

    def get_success_url(self):
        messages.success(self.request, "Gönderi başarılı bir şekilde silindi")
        return reverse('blog:anasayfa')


class CategoryDeleteView(generic.DeleteView):
    template_name = "blog/pages/delete.html"
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
    template_name = "blog/pages/comments.html"
    model = Comments

    def get(self, request, *args, **kwargs):
        if not request.user.is_staff:
            messages.error(request, "Yetkili girişi yapınız!")
            return redirect('%s?next=/blog/' % (settings.LOGIN_URL))
        return super().get(request, *args, **kwargs)


def comment_read(request, pk):
    if request.user.is_anonymous:
        messages.error(request, "Gönderiyi beğenmek için giriş yapınız!")
        return redirect('%s?next=/blog/' % (settings.LOGIN_URL))
    if not request.user.is_staff:
        messages.error(request, "Yetkili girişi yapınız!")
        return redirect('%s?next=/blog/' % (settings.LOGIN_URL))
    comment = get_object_or_404(Comments, pk=pk)
    comment.status = "okundu"
    comment.save()
    return HttpResponseRedirect(reverse("dashboard:anasayfa", args=[request.user.username]))


def comment_delete(request, pk):
    if request.user.is_anonymous:
        messages.error(request, "Gönderiyi beğenmek için giriş yapınız!")
        return redirect('%s?next=/blog/' % (settings.LOGIN_URL))
    comment = get_object_or_404(Comments, pk=pk)
    if not request.user.username == comment.commentator.username:
        messages.error(request, "Yetkisiz kullanıcı")
        return redirect('%s?next=/blog/' % (settings.LOGIN_URL))
    comment.delete()
    messages.success(request, "Yorumunuz başarılı bir şekilde silindi")
    return HttpResponseRedirect(request.META.get('HTTP_REFERER'))


def comment_like(request, id):
    if request.user.is_anonymous:
        messages.error(request, "Gönderiyi beğenmek için giriş yapınız!")
        return redirect('%s?next=/blog/' % (settings.LOGIN_URL))
    comment = get_object_or_404(Comments, id=id)

    if comment.commentator.username == request.user.username:
        messages.info(request, "Kendini beğenmiş birisin :) ")
    else:
        messages.success(request, f"{comment.commentator.username} isimli kullanıcının yorumunu beğendiniz")

    comment.likes.add(request.user)

    if not Notifications.objects.filter(user=comment.commentator, title=f"{request.user} yorumunu beğendi",
                                 type="dislike_comment", text=comment.comment):

        Notifications.objects.create(user=comment.commentator, title=f"{request.user} yorumunu beğendi",
                                     type="dislike_comment", text=comment.comment)

        for like in comment.likes.values():
            if not request.user.username == like.get("username"):
                send_mail(
                    f"{request.user} yorumunuzu beğendi",
                    f"{request.user} isimli kullanıcı {comment.comment} yorumunuzu beğendi.",
                    "bioblogdestek@gmail.com",
                    [comment.commentator.email],
                    fail_silently=False,
                )

    return HttpResponseRedirect(request.META.get('HTTP_REFERER'))


def comment_dislike(request, id):
    if request.user.is_anonymous:
        messages.error(request, "Gönderiyi beğenmek için giriş yapınız!")
        return redirect('%s?next=/blog/' % (settings.LOGIN_URL))

    comment = get_object_or_404(Comments, id=id)
    comment.dislike.add(request.user)
    messages.error(request, f"{comment.commentator.username} isimli kullanıcının yorumunu beğenmediniz")

    if not Notifications.objects.filter(
            user=comment.commentator,
            title=f"{request.user} yorumunu beğenmedi",
            type="dislike_comment", text=comment.comment):

        Notifications.objects.create(user=comment.commentator, title=f"{request.user} yorumunu beğenmedi",
                                     type="dislike_comment", text=comment.comment)

        for like in comment.dislike.values():
            if not request.user.username == like.get("username"):

                send_mail(
                    f"{request.user} yorumunuzu beğenmedi",
                    f"{request.user} isimli kullanıcı {comment.comment} yorumunuzu beğenmedi.",
                    "bioblogdestek@gmail.com",
                    [comment.commentator.email],
                    fail_silently=False,
                )

    return HttpResponseRedirect(request.META.get('HTTP_REFERER'))


def report_comment(request, id):
    if request.user.is_anonymous:
        messages.error(request, "Gönderiyi beğenmek için giriş yapınız!")
        return redirect('%s?next=/blog/' % (settings.LOGIN_URL))
    comment = get_object_or_404(Comments, id=id)
    comment.report.add(request.user)
    if comment.report.count() > 5:
        comment.delete()
    messages.success(request, "Geri bildiriminiz iletilmiştir")
    return HttpResponseRedirect(request.META.get('HTTP_REFERER'))


def like_post(request, pk, slug):
    if request.user.is_anonymous:
        messages.error(request, "Gönderiyi beğenmek için giriş yapınız!")
        return redirect('%s?next=/blog/' % (settings.LOGIN_URL))

    post = get_object_or_404(Posts, pk=pk, slug=slug)
    post.likes.add(request.user)
    messages.success(request, f"{post.title} başlıklı gönderiyi beğendiniz")

    if not Notifications.objects.filter(user=post.author, title=f"{request.user}, {post.title} başlıklı gönderinizi beğendi.",
                                 type="like_post"):
        Notifications.objects.create(user=post.author,
                                     title=f"{request.user}, {post.title} başlıklı gönderinizi beğendi.",
                                     type="like_post")
    for like in post.likes.values():
        if not request.user.username == like.get("username"):
            send_mail(
                f"{request.user} gönderinizi beğendi",
                f"{request.user} isimli kullanıcı {post.title} başlıklı gönderinizi beğendi.",
                "bioblogdestek@gmail.com",
                [post.author.email],
                fail_silently=False,
            )
    return HttpResponseRedirect(request.META.get('HTTP_REFERER'))


def dislike_post(request, pk, title):
    if request.user.is_anonymous:
        messages.error(request, "Gönderiyi beğenmek için giriş yapınız!")
        return redirect('%s?next=/blog/' % (settings.LOGIN_URL))
    post = get_object_or_404(Posts, pk=pk, title=title)
    post.dislike.add(request.user)
    messages.error(request, f"{post.title} başlıklı gönderiyi beğenmediniz")
    if not Notifications.objects.filter(user=post.author,
                                 title=f"{request.user}, {post.title} başlıklı gönderinizi beğenmedi.",
                                 type="dislike_post"):
        Notifications.objects.create(user=post.author,
                                     title=f"{request.user}, {post.title} başlıklı gönderinizi beğenmedi.",
                                     type="dislike_post")

    return HttpResponseRedirect(request.META.get('HTTP_REFERER'))


class CommentUpdate(generic.UpdateView):
    template_name = "blog/pages/update.html"
    model = Comments
    form_class = AddCommentForm

    def get_success_url(self):
        obj = self.get_object().post.pk
        post = get_object_or_404(Posts, pk=int(obj))
        print(post)
        messages.success(self.request, "Yorumunuz başarılı bir şekilde güncellendi..")
        return reverse('blog:post_detail',
                       kwargs={
                           'category': post.category.slug,
                           'slug': post.slug,
                           'pk': post.pk,
                           'author': post.author,
                           'author_id': post.author_id,
                           'created': post.created,
                       })


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
        messages.success(request, "Mesajınız yazara iletildi")
        return HttpResponseRedirect(self.request.build_absolute_uri())

    def get_context_data(self, *, object_list=None, **kwargs):
        context = super().get_context_data(**kwargs)
        for profile in Profile.objects.filter(user=self.kwargs['pk']):
            context['user_social'] = profile.user_social.values()
            context['posts_lists'] = Posts.objects.filter(author=self.kwargs['pk'])
        return context


def profile_view(request, username, pk):
    global user_social, object
    post_list = Posts.objects.filter(author=pk)

    try:
        object = User.objects.get(username=username, pk=pk)

    except:
        return render(request, "exception/page-404.html")

    for profile in Profile.objects.filter(user=pk):
        user_social = profile.user_social.values()

    form = ContactProfileForm(request.POST or None)

    # Paginator
    paginator = Paginator(post_list, 6)
    page_number = request.GET.get("page")
    page_obj = paginator.get_page(page_number)

    if request.method == "POST":

        email = request.POST.get('email')

        if email:

            if Subscribe.objects.exists():
                for abone in Subscribe.objects.all():
                    if email in abone.email:
                        messages.error(request, "Bültene zaten abonesiniz")
                        return HttpResponseRedirect(request.build_absolute_uri())

                    else:
                        Subscribe.objects.create(email=email)

            else:
                Subscribe.objects.create(email=email)

        if form.is_valid():

            content = form.cleaned_data['content']
            contact_email = form.cleaned_data['contact_email']
            title = form.cleaned_data['title']
            receiver = User.objects.get(pk=pk)
            sender = request.user

            if sender.username == receiver.username:
                messages.error(request, "Kendi kendinize mesaj attınız")

                return HttpResponseRedirect(request.build_absolute_uri())

            ContactModel.objects.create(sender=sender, receiver=receiver, content=content, contact_email=contact_email,
                                        title=title)

            send_mail(
                "Mesajınız var.",
                f" Gönderen : {sender} \n\n\n Başlık: {title} \n\n\n ,"
                f"Email: {contact_email} \n\n\n Yazdığı mesaj : \n\n\n {content}",
                "bioblogdestek@gmail.com",
                [receiver.email],
                fail_silently=False,
            )

            messages.success(request, "Mesajınız başarılı bir şekilde gönderildi")

            return HttpResponseRedirect(request.build_absolute_uri())

    return render(request, 'blog/pages/profile.html', {'form': form, 'post_list': post_list,
                                                       'user_social': user_social, "object": object,
                                                       "page_obj": page_obj})


class ProfileUpdateViewNonStaff(generic.UpdateView):
    template_name = 'blog/pages/update.html'
    model = Profile
    form_class = ProfileModelForm

    def form_valid(self, form):
        instance = form.save(commit=False)
        instance.user = self.request.user
        instance.save()
        messages.success(self.request, "Profil güncellendi!")
        return HttpResponseRedirect(reverse('blog:profile', kwargs={
            'username': self.request.user, 'pk': self.request.user.pk
        }))


def blog_contact(request):
    contact = BlogContactModel.objects.last()
    contact_form = BlogContactForm(request.POST or None)

    if request.method == "POST":
        email = request.POST.get('email_sub')

        if email:
            if Subscribe.objects.exists():
                for abone in Subscribe.objects.all():
                    if email in abone.email:
                        messages.error(request, "Bültene zaten abonesiniz")
                        return HttpResponseRedirect(request.build_absolute_uri())

                    else:
                        Subscribe.objects.create(email=email)
                        messages.success(request, "Başarılı bir şekilde abone oldunuz")

                        return HttpResponseRedirect(request.build_absolute_uri())

            else:
                Subscribe.objects.create(email=email)

                messages.success(request, "Başarılı bir şekilde abone oldunuz")

                return HttpResponseRedirect(request.build_absolute_uri())

        elif contact_form.is_valid():
            name = contact_form.cleaned_data["name"]
            email = contact_form.cleaned_data["email"]
            message = contact_form.cleaned_data["message"]

            send_mail(
                "Blog iletişime mesaj atan oldu...",
                f" Gönderen : {name} \n\n\n Email: {email} \n\n\n Yazdığı mesaj : \n\n {message}",
                "bioblogdestek@gmail.com",
                ["ozer246@gmail.com"],
                fail_silently=False,
            )

            messages.success(request, "Mesajınız iletildi. En kısa sürede mesajınızı yanıtlamaya çalışacağız. ")

            contact_form.save()

            return HttpResponseRedirect(request.build_absolute_uri())

        else:
            contact_form = BlogContactForm()

    return render(request, "blog/pages/blog_contact.html", {'contact': contact, 'contact_form': contact_form})


def about(request):
    if request.method == "POST":
        email = request.POST.get('email_sub')

        if email:

            if Subscribe.objects.exists():

                for abone in Subscribe.objects.all():

                    if email in abone.email:
                        messages.error(request, "Bültene zaten abonesiniz")
                        return HttpResponseRedirect(request.build_absolute_uri())

                    else:
                        Subscribe.objects.create(email=email)
                        messages.success(request, "Başarılı bir şekilde abone oldunuz")

                        return HttpResponseRedirect(request.build_absolute_uri())

            else:

                Subscribe.objects.create(email=email)

                messages.success(request, "Başarılı bir şekilde abone oldunuz")

                return HttpResponseRedirect(request.build_absolute_uri())

    return render(request, "blog/pages/about.html")


def abone_ol(request):
    email = request.POST.get('email_sub')
    if Subscribe.objects.exists():
        for abone in Subscribe.objects.all():
            if email in abone.email:
                messages.error(request, "Bültene zaten abonesiniz")
                return redirect(request.META['HTTP_REFERER'])
            else:
                Subscribe.objects.create(email=email)
                messages.success(request, "Başarılı bir şekilde abone oldunuz")
                return redirect(request.META['HTTP_REFERER'])
    else:
        Subscribe.objects.create(email=email)
        messages.success(request, "Başarılı bir şekilde abone oldunuz")
        return redirect(request.META['HTTP_REFERER'])
